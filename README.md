# My Technical Blog

A personal technical blog built with Jekyll and the Chirpy theme, hosted on GitHub Pages.

## About

This blog shares insights on software development, programming, and technology. Topics covered include:

- Web Development
- Programming Languages (Python, JavaScript, etc.)
- Version Control (Git)
- DevOps & Best Practices
- And more!

## Live Site

Visit the blog at: [https://myx0506.github.io](https://myx0506.github.io)

## Local Development

To run this blog locally:

1. Install Ruby and Bundler
2. Clone this repository
3. Install dependencies:
   ```bash
   bundle install
   ```
4. Run the development server:
   ```bash
   bundle exec jekyll serve
   ```
5. Open your browser to `http://localhost:4000`

## Writing Posts

Create new posts in the `_posts` directory following this naming convention:
```
YYYY-MM-DD-title-of-post.md
```

Each post should have front matter at the top:
```yaml
---
title: Your Post Title
date: YYYY-MM-DD HH:MM:SS -0500
categories: [Category1, Category2]
tags: [tag1, tag2, tag3]
---
```

## Customization

- Edit `_config.yml` to update site settings
- Modify `_tabs/about.md` to update the About page
- Add your own avatar image and update the `avatar` field in `_config.yml`

## Deployment

This site automatically deploys to GitHub Pages when changes are pushed to the `main` or `master` branch.

## Theme

This blog uses the [Chirpy](https://github.com/cotes2020/jekyll-theme-chirpy) theme, a minimal, responsive, and feature-rich Jekyll theme for technical writing.

## License

This work is published under [MIT](LICENSE) License.
